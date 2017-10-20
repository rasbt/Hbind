#include <string.h>
#include "analyze_ligand.h"
#include "adj_list.h"
#include "find_carbon_ring_centers.h"
#include "assign_hydrogens.h"
#include "find_cycles.h"
#include "find_flexible_bonds.h"
#include "assign_fragments.h"
#include "number_ligand_atoms.h"

int  analyze_ligand ( global_data_pt  global )
{
  pts_compound_pt compounds;
  molecule_pt     ligand;
  int             *atom_index,
                  *act,
                  *neighbors;
  int             first,
                  real_index;
  int             i, j;

  ligand = global->ligand;
  compounds = global->compound_interactions->compounds;
  atom_index = global->compound_interactions->atom_index;
  act = global->compound_interactions->act;
 
  construct_adjacency_list(ligand);
  find_carbon_ring_centers(ligand);
  if(assign_hydrogens(ligand) == FAILURE) return FAILURE;  
  if(find_cycles(ligand) == FAILURE) return FAILURE;
  if(find_flexible_bonds(ligand, global->flex_bond_rules,
			 global->number_of_flex_bond_rules) == FAILURE) 
    return FAILURE;
  number_ligand_atoms(ligand);
  assign_fragments(ligand);

  first = compounds[global->compound_index].first_point;
  for ( i = 0; i < compounds[global->compound_index].number_of_points; i++ ){

    if ( act[first+i] == HYDROPHOB ) continue;
      
    real_index = ligand->atom_index[atom_index[first+i]+1];
    ligand->atoms[real_index].act = act[first+i];
    if ( act[first+i] == DONOR || act[first+i] == DONEPTOR ){
      neighbors = ligand->neighbors[real_index];
      for ( j = 0; j < ligand->number_of_neighbors[real_index]; j++ )
        if ( ligand->atoms[neighbors[j]].type == H )
	  ligand->atoms[neighbors[j]].act = DONOR_HYDROGEN;
    }
  }
  memcpy(ligand->orig_atom_positions, ligand->atom_positions, 
         3*ligand->number_of_atoms * sizeof(*ligand->atom_positions));
  for ( i = 0; i < ligand->number_of_atoms; i++ ) {
    ligand->orig_act[i] = ligand->atoms[i].act;
    ligand->orig_rad[i] = ligand->atoms[i].rad;
  }
  return SUCCESS;
}

