#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "bitstrings.h"
#include "mymalloc.h"

global_data_pt  initialize_global_data_structure ( void )
{
  global_data_pt  global;
  int             i;

  global = (global_data_pt) mymalloc ( sizeof (global_data_t) );

  global->ligand_flag = 0;
  global->match_2_key_points = FALSE;
  global->total_num_molecules_noconf = 0;
  global->total_num_molecules_conf = 0;
  global->total_num_output_molecules = 0;
  global->best_affiscore_so_far = 0.0;
  
  global->hbind_dir = 
    (char *) mymalloc ( FILENAME_MAX );
  global->data_root = 
    (char *) mymalloc ( FILENAME_MAX );
  *global->hbind_dir = '\0';
  *global->data_root = '\0';
  global->compound_name = 
    (char *) mymalloc ( FILENAME_MAX );
  global->compound_dir = 
    (char *) mymalloc ( FILENAME_MAX );
  global->protein = 
    (char *) mymalloc ( FILENAME_MAX );
  global->template = 
    (char *) mymalloc ( FILENAME_MAX );
  global->database = 
    (char *) mymalloc ( FILENAME_MAX );
  global->best_affi_name =    
    (char *) mymalloc ( FILENAME_MAX );
    
  /* allocate vectors for data of target protein */
  global->target_atoms = 
    (atom_pt) mymalloc ( MAX_PDB_ATOMS * sizeof (atom_t) );
  global->target_residues = 
    (residue_pt) mymalloc ( MAX_PDB_RESIDUES * sizeof (residue_t) );
  global->target_rotations = 
    (int *) mymalloc ( MAX_PDB_RESIDUES * sizeof (int) );
  global->target_intra_overlap = 
    (int *) mymalloc ( MAX_PDB_RESIDUES * sizeof (int) );
  global->target_atom_positions = 
    (float *) mymalloc(3*MAX_PDB_ATOMS*sizeof(*global->target_atom_positions));
  global->orig_target_atom_positions =
    (float *) mymalloc(3*MAX_PDB_ATOMS * sizeof(float));
  global->orig_target_atom_act =
    (int *) mymalloc ( MAX_PDB_ATOMS * sizeof (int) );
  for(i = 0; i < MAX_PDB_ATOMS; i++)
    global->target_atoms[i].pos = &(global->target_atom_positions[3*i]);

  global->template_interactions =
    (interaction_pt) mymalloc ( MAX_TEMPLATE_POINTS * sizeof (interaction_t));
 
  global->compound_interactions = (pts_pt) mymalloc ( sizeof (pts_t) );
  global->compound_interactions->compounds =
    (pts_compound_pt) mymalloc ( MAX_NUMBER_OF_PTS_COMPOUNDS 
				 * sizeof (pts_compound_t) );
  global->compound_interactions->pos =
    (float *) mymalloc ( 3 
			 * MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS
			 * sizeof (float) );
  global->compound_interactions->act =
    (int *) mymalloc ( MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS 
		       * sizeof (int) );
  global->compound_interactions->atom_index =
    (int *) mymalloc ( MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS 
		       * sizeof (int) );
  
  global->ligand_file_name = 
    (char *) mymalloc ( FILENAME_MAX * sizeof (char) );
  global->old_compound_name = 
    (char *) mymalloc ( FILENAME_MAX * sizeof (char) );
  global->old_compound_name[0] = '\0';

  global->old_ligand_name_noconf = 
    (char *) mymalloc ( FILENAME_MAX * sizeof (char) );
  global->old_ligand_name_noconf[0] = '\0';

  global->flex_bond_rules = 
    (flex_bond_defn_t *) mymalloc ( MAX_NUMBER_OF_FLEX_BOND_RULES 
				    * sizeof (flex_bond_defn_t) );


  global->ligand_bumps = 
    (int *) mymalloc ( MAX_SIDE_CHAIN_BUMPS * sizeof (int) );
  global->target_bumps = 
    (int *) mymalloc ( MAX_SIDE_CHAIN_BUMPS * sizeof (int) );
  global->unbump_target_indices =
    (int **) mymalloc ( MAX_PDB_RESIDUES * sizeof (int *) );
  for ( i = 0; i < MAX_PDB_RESIDUES; i++ )
    global->unbump_target_indices[i] = 
      (int *) mymalloc ( 5 * sizeof (int) );
  global->unbump_indices = 
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  global->unbump_bonds = 
    (int *) mymalloc ( MAX_UNBUMP_BONDS * sizeof (int) );
  global->unbump_matrices =
    (transform_matrix_pt) mymalloc ( MAX_SIDE_CHAIN_BUMPS * MAX_UNBUMP_BONDS
				     * sizeof (transform_matrix_t) );
  global->unbump_data = 
    (unbump_data_pt **) mymalloc ( MAX_SIDE_CHAIN_BUMPS
				  * sizeof (unbump_data_pt *) );
  for ( i = 0; i < MAX_SIDE_CHAIN_BUMPS; i++ )
    global->unbump_data[i] = 
      (unbump_data_pt *) mymalloc ( MAX_UNBUMP_BONDS 
				  * sizeof (unbump_data_pt) );
  global->unbump_entries = 
    (unbump_data_pt) mymalloc ( MAX_SIDE_CHAIN_BUMPS * MAX_UNBUMP_BONDS / 3
				* sizeof (unbump_data_t) );
  global->unbump_dependencies =
    (int *) mymalloc ( MAX_UNBUMP_BONDS 
		       * MAX_UNBUMP_DEPENDENCIES 
		       * sizeof (int) );
  global->unbump_dependents = 
    (int **) mymalloc ( MAX_UNBUMP_BONDS * sizeof (int *) );
  global->number_of_unbump_dependents = 
    (int *) mymalloc ( MAX_UNBUMP_BONDS * sizeof (int) );

  global->ligand = (molecule_pt) mymalloc ( sizeof (molecule_t) );
  global->ligand->name[0] = '\0';
  global->ligand->name_noconf[0] = '\0';
  global->ligand->conf_number[0] = '\0';
  global->ligand->atoms = 
    (atom_pt) mymalloc(MAX_NUMBER_OF_MOL2_ATOMS * sizeof(atom_t));
  global->ligand->atom_positions = 
    (float*) mymalloc(3*MAX_NUMBER_OF_MOL2_ATOMS * sizeof(float));
  global->ligand->orig_atom_positions = 
    (float *) mymalloc (3*MAX_NUMBER_OF_MOL2_ATOMS * sizeof(float));
  global->ligand->orig_act = 
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );
  global->ligand->orig_rad = 
    (float *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (float) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    global->ligand->atoms[i].pos = &(global->ligand->atom_positions[3*i]);

  /** Added by RSK - for substructures - Aug 13 2002 **/
  for ( i = 0; i < MAX_SUBSTS; i++ )
  global->ligand->substructure[i] = 
      (char *) mymalloc(MAX_SUBST_LEN * sizeof(char ));

  global->ligand->relations = 
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS
		       * MAX_NUMBER_OF_MOL2_ATOMS
		       * sizeof (int) );
  global->ligand->atom_index =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );
  global->ligand->bonds = 
    (bond_pt) mymalloc ( MAX_NUMBER_OF_MOL2_BONDS * sizeof (bond_t) );
  global->ligand->flexible_bonds = 
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  global->ligand->bond_types = 
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  global->ligand->bond_directions = 
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  global->ligand->fragment_locations = 
    (int *) mymalloc ( ( MAX_NUMBER_OF_FLEXIBLE_BONDS + 1 ) * sizeof (int) );
  global->ligand->anchor_dist = 
    (int *) mymalloc ( ( MAX_NUMBER_OF_FLEXIBLE_BONDS + 1 ) * sizeof (int) );
  global->ligand->way_to_anchor = 
    (int *) mymalloc ( ( MAX_NUMBER_OF_FLEXIBLE_BONDS + 1 ) * sizeof (int) );
  global->ligand->carbon_ring_centers = 
    (float **) mymalloc ( MAX_NUMBER_OF_CARBON_RINGS * sizeof (float *) );
  for ( i = 0; i < MAX_NUMBER_OF_CARBON_RINGS; i++ )
    global->ligand->carbon_ring_centers[i] =
      (float *) mymalloc ( 3 * sizeof (float) );
  global->ligand->carbon_ring_atom = 
    (int *) mymalloc ( MAX_NUMBER_OF_CARBON_RINGS * sizeof (int) );
  global->ligand->neighbors = 
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    {
      global->ligand->neighbors[i] =
	(int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
      global->ligand->atoms[i].fragments = 
	bitstring_create ( MAX_NUMBER_OF_FLEXIBLE_BONDS );
      bitstring_clear_all ( global->ligand->atoms[i].fragments );
    }
  global->ligand->bond_to_neighbor = 
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    global->ligand->bond_to_neighbor[i] =
      (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
  global->ligand->number_of_neighbors =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );
  global->ligand->bond_order =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );

  global->ligand->fragment_neighbors = 
    (int **) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_FLEXIBLE_BONDS; i++ )
    global->ligand->fragment_neighbors[i] =
      (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );

  global->ligand->bond_to_fragment_neighbor = 
    (int **) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_FLEXIBLE_BONDS; i++ )
    global->ligand->bond_to_fragment_neighbor[i] =
      (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  global->ligand->number_of_fragment_neighbors =
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );

  /* Get triangle parameters from defs.h */
  global->triangle_parameters.min_perimeter = TRIANGLE_MIN_PERIMETER;
  global->triangle_parameters.max_perimeter = TRIANGLE_MAX_PERIMETER;
  global->triangle_parameters.min_short_side = TRIANGLE_MIN_SHORTEST_SIDE;
  global->triangle_parameters.max_short_side = TRIANGLE_MAX_SHORTEST_SIDE;
  global->triangle_parameters.max_long_side = TRIANGLE_MAX_LONGEST_SIDE;
  /* This is set to the smaller of the 2 values because for one screening run
   * we may have large and small ligands.  This means that the hash table must
   * include buckets corresponding to the union of all the ranges.  */
  global->triangle_parameters.min_long_side = SMALL_TRIANGLE_MIN_LONGEST_SIDE;

  global->MM2_FILE = NULL; 

  global->number_of_screened_compounds = 0;
  global->number_of_potential_ligands = 0;
  for ( i = 0; i < NUMBER_OF_FILTERS; i++ )
    global->filter_counter[i] = 0;
  global->number_of_bumped_ligands = 0;
  global->score_only = FALSE;
  global->restart_molecule_check = FALSE;
  global->restart_molecule[0] = 0;
  return global;
}
