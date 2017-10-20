#include <math.h>
#include <stdio.h>
#include <mymalloc.h>
#include "types.h"
#include "dist_fun.h"
#include "bitstrings.h"
#include "unbump_anchor.h"
#include "unbump_side_chains.h"
#include "intra_bump_check.h"
#include "find_all_bumps.h"
#include "calc_score_from_terms.h"
#include "docking_features.h"
#include <distance_matrices.h>

#define AFFINITY_SCORE 62  		/* affinity trained funtion # 66 */
#define ORIENTATION_SCORE 63	/* affinity trained funtion # 66 */


void ligand_trace(molecule_pt ligand);

void print_target_ligand_distances(atom_pt target, int num_target_atoms,
                                   atom_pt ligand, int num_ligand_atoms);

void traverse_fragment_graph(molecule_pt ligand, int  fragment, int neighbor)
{
  int old_fragment;
  int i;
  int neighbor_bond = ligand->bond_to_fragment_neighbor[fragment][neighbor];

  /* this bond is located inside the anchor fragment */
  if(ligand->bond_types[neighbor_bond] == RIGID )
    ligand->fragment_locations[fragment] = ANCHOR;
  old_fragment = fragment;
  fragment = ligand->fragment_neighbors[fragment][neighbor];
  for(i = 0; i < ligand->number_of_fragment_neighbors[fragment]; i++ ){
    neighbor = ligand->fragment_neighbors[fragment][i];
    if ( neighbor != old_fragment
	 && ligand->fragment_locations[neighbor] != ANCHOR )
      traverse_fragment_graph ( ligand, fragment, i );
  }
}

/* We already know the three (or less) fragments that contain the mapped
 * ligand interaction centers global->anchor_fragments, and we know those
 * flexible bonds that were rigidified, since they connect these fragments.
 * The goal of this routine is to mark the full anchor fragment and the
 * anchor atoms, these are those atoms that cannot be displaced, since
 * they are directly connected to the anchor fragment by a flexible bond.
 */
void  find_anchor_fragments_and_anchor_atoms ( global_data_pt  global )
{
  molecule_pt ligand;
  bond_pt     bond;
  int         fragment;
  int         i, j;
  
  ligand = global->ligand;
  for ( i = 0; i < ligand->number_of_fragments; i++ )
    ligand->fragment_locations[i] = OFF_ANCHOR;
  for ( i = 0; i < 3; i++ )
    ligand->fragment_locations[global->anchor_fragments[i]] = ANCHOR;
  for ( i = 0; i < 3; i++ ) {
    fragment = global->anchor_fragments[i];
    for ( j = 0; j < ligand->number_of_fragment_neighbors[fragment]; j++ )
      traverse_fragment_graph ( ligand, fragment, j );
  }
  /* so far, we set the array fragment_location[] to the correct values, now
     the atoms have to be marked */
  for ( i = 0; i < ligand->number_of_atoms; i++ )
    if ( ligand->fragment_locations[ligand->atoms[i].fragment] == ANCHOR )
      ligand->atoms[i].res = ANCHOR;
    else ligand->atoms[i].res = OFF_ANCHOR;

  /* we have to mark some additional atoms as anchor-fragment atoms, these
     are those atoms that are adjacent to a flexible bond that is connected
     to the anchor fragment (similar to C-beta atoms in proteins, which
     technically are side-chain atoms, but with a rigid main chain, there
     position is fixed) */
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( ligand->bond_types[i] == FLEX ) {
      bond = &ligand->bonds[ligand->flexible_bonds[i]];
      if ( ligand->atoms[bond->atom1].res == ANCHOR )
	ligand->atoms[bond->atom2].res = ANCHOR;
      if ( ligand->atoms[bond->atom2].res == ANCHOR )
	ligand->atoms[bond->atom1].res = ANCHOR;
    }
}

int check_buried_carbons(global_data_pt global)
{
  int i;
  atom_pt ligand = global->ligand->atoms;
  int number_of_carbons = 0;
  int number_of_contact_carbons = 0;
  dist_bin_pt bin = 0;
  atom_pt *targ_atom_p;

  for(i = 0; i < global->ligand->number_of_atoms; i++ ){
    if(ligand[i].type == C ) number_of_carbons++;
    bin = get_bin(&global->target_dists_array, ligand[i].pos);
    if(bin == 0) continue;
    for(targ_atom_p = bin->atoms; targ_atom_p < bin->atoms + bin->num_atoms; 
        ++targ_atom_p){
      if(squared_dist(ligand[i].pos, (*targ_atom_p)->pos) < HYDRO_DIST_2){
        number_of_contact_carbons++;
        break;
      }
    }
  }
  if(number_of_carbons > 0 
     && (float) number_of_contact_carbons / (float) number_of_carbons > 0.5 )
    return SUCCESS;
  return FAILURE;
}


void find_rigid_bonds(molecule_pt ligand, int *ligand_matches)
{
  int i;
  bitstring_pt strings[3];
  atom_pt atoms = ligand->atoms;
  
  for ( i = 0; i < 3; i++ )
    strings[i] = bitstring_create ( MAX_NUMBER_OF_FLEXIBLE_BONDS );
  bitstring_cpy_all ( strings[0], atoms[ligand_matches[0]].fragments );
  bitstring_xor_all ( strings[0], atoms[ligand_matches[1]].fragments );
  bitstring_cpy_all ( strings[1], atoms[ligand_matches[1]].fragments );
  bitstring_xor_all ( strings[1], atoms[ligand_matches[2]].fragments );
  bitstring_cpy_all ( strings[2], atoms[ligand_matches[2]].fragments );
  bitstring_xor_all ( strings[2], atoms[ligand_matches[0]].fragments );
  bitstring_or_all ( strings[0], strings[1] );
  bitstring_or_all ( strings[0], strings[2] );
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( bitstring_get_bit ( strings[0], i ) ) ligand->bond_types[i] = RIGID;
    else ligand->bond_types[i] = FLEX;
  for ( i = 0; i < 3; i++ ) bitstring_free ( strings[i] );
}  

void  direct_bond ( molecule_pt  molecule,
		    int          *mark,
		    int          flex_bond,
		    int          dist,
		    int          start_atom )
{
  atom_pt  atoms;
  bond_pt  bonds;
  int      **fragment_bonds;
  int      *flex_bonds;
  int      fragment,
    number_of_neighbors,
    next_bond,
    bond_index;
  int      i;
  
  atoms = molecule->atoms;
  bonds = molecule->bonds;
  flex_bonds = molecule->flexible_bonds;
  bond_index = flex_bonds[flex_bond];
  fragment_bonds = molecule->bond_to_fragment_neighbor;
  if ( start_atom == ATOM1 ) {
    fragment = atoms[bonds[bond_index].atom2].fragment;
    mark[fragment] = VISITED;
    molecule->anchor_dist[fragment] = dist;
    molecule->way_to_anchor[fragment] = flex_bond;
    number_of_neighbors = molecule->number_of_fragment_neighbors[fragment];
    for ( i = 0; i < number_of_neighbors; i++ )
      /* don't go in cycles */
      if ( mark[molecule->fragment_neighbors[fragment][i]] == UNVISITED ){
	next_bond = fragment_bonds[fragment][i];
	if ( atoms[bonds[flex_bonds[next_bond]].atom1].fragment == fragment )
	  direct_bond ( molecule, mark, next_bond, dist + 1, ATOM1 );
	else direct_bond ( molecule, mark, next_bond, dist + 1, ATOM2 );
      }
    /* numbering hadn't changed order of adjacent atoms */
    if ( atoms[bonds[bond_index].atom1].number 
	 < atoms[bonds[bond_index].atom2].number )
      molecule->bond_directions[flex_bond] = REVERSE;
    /* numbering had changed order, so direction is opposed to
       what one expects when only checking .atom1 and .atom2 */
    else molecule->bond_directions[flex_bond] = STRAIGHT;
  } else {
    fragment = atoms[bonds[bond_index].atom1].fragment;
    mark[fragment] = VISITED;
    molecule->anchor_dist[fragment] = dist;
    molecule->way_to_anchor[fragment] = flex_bond;
    number_of_neighbors = molecule->number_of_fragment_neighbors[fragment];
    for ( i = 0; i < number_of_neighbors; i++ )
      /* don't go in cycles */
      if ( mark[molecule->fragment_neighbors[fragment][i]] == UNVISITED ) {
	next_bond = fragment_bonds[fragment][i];
	if ( atoms[bonds[flex_bonds[next_bond]].atom1].fragment == fragment )
	  direct_bond ( molecule, mark, next_bond, dist + 1, ATOM1 );
	else direct_bond ( molecule, mark, next_bond, dist + 1, ATOM2 );
      }
    /* numbering hadn't changed order of adjacent atoms */
    if ( atoms[bonds[bond_index].atom1].number 
	 < atoms[bonds[bond_index].atom2].number )
      molecule->bond_directions[flex_bond] = STRAIGHT;
    else
      molecule->bond_directions[flex_bond] = REVERSE;
  }
}

void  direct_flexible_bonds ( global_data_pt  global )
{
  molecule_pt molecule;
  bond_pt     bonds;
  atom_pt     atoms;
  int         **neighbors,
    **bond;
  int         *flex_bonds,
    *number_of_neighbors;
  int         mark[MAX_NUMBER_OF_FLEXIBLE_BONDS];
  int         fragment1,
    fragment2;
  int         i;
  
  molecule = global->ligand;
  atoms = molecule->atoms;
  flex_bonds = molecule->flexible_bonds;
  bonds = molecule->bonds;
  neighbors = molecule->fragment_neighbors;
  number_of_neighbors = molecule->number_of_fragment_neighbors;
  for ( i = 0; i < molecule->number_of_fragments; i++ ) {
    number_of_neighbors[i] = 0;
    molecule->anchor_dist[i] = 0;
    molecule->way_to_anchor[i] = UNKNOWN;
  }
  bond = molecule->bond_to_fragment_neighbor;
  for ( i = 0; i < molecule->number_of_flexible_bonds; i++ ) {
    fragment1 = atoms[bonds[flex_bonds[i]].atom1].fragment;
    fragment2 = atoms[bonds[flex_bonds[i]].atom2].fragment;
    neighbors[fragment1][number_of_neighbors[fragment1]] = fragment2;
    bond[fragment1][number_of_neighbors[fragment1]++] = i;
    neighbors[fragment2][number_of_neighbors[fragment2]] = fragment1;
    bond[fragment2][number_of_neighbors[fragment2]++] = i;
  }
  
#ifdef TRACE
  int j;
  for ( i = 0; i < molecule->number_of_fragments; i++ ) {
    printf ( "Fragment %d: ", i );
    for ( j = 0; j < number_of_neighbors[i]; j++ )
      printf ( "%d (%d) ", neighbors[i][j] + 1, flex_bonds[bond[i][j]] );
    printf ( "\n" );
  }
#endif
  
  for ( i = 0; i < molecule->number_of_fragments; i++ ) mark[i] = UNVISITED;
  for ( i = 0; i < molecule->number_of_flexible_bonds; i++ )
    if ( molecule->bond_types[i] == FLEX ) {
      if(molecule->fragment_locations[atoms[bonds[flex_bonds[i]].atom1].
				      fragment] == ANCHOR ) {
	mark[atoms[bonds[flex_bonds[i]].atom1].fragment] = VISITED;
	molecule->anchor_dist[atoms[bonds[flex_bonds[i]].atom1].fragment] 
	  = 0;
	direct_bond ( global->ligand, mark, i, 1, ATOM1 );
      }else if(molecule->fragment_locations[atoms[bonds[flex_bonds[i]].atom2].
					    fragment] == ANCHOR ) {
	mark[atoms[bonds[flex_bonds[i]].atom2].fragment] = VISITED;
	molecule->anchor_dist[atoms[bonds[flex_bonds[i]].atom2].fragment] = 0;
	direct_bond ( global->ligand, mark, i, 1, ATOM2 );  
      }
    }
}

/*
 * This routine is the last one that is called before scoring the complex.
 * It's main purpose is to check, if any bumps slipped through, which shouldn't
 * happen. In addition to this, all waters are checked for collisions with any
 * ligand or protein atom. There might be some collsions, because waters are
 * not considered during intramolecular bump checks. 
 */
int full_bump_check ( global_data_pt global )
{
  atom_pt  water;
  float distance;
  float overlap;
  float sq_dist;
  int i, j;
  float max_overlap = global->finally_tolerated_max_bump;  
  int resolve = SUCCESS;
  int counter = 0;
  const int num_targ_atoms = global->number_of_target_atoms;
  atom_pt ligand = global->ligand->atoms;
  atom_pt target = global->target_atoms;

  for ( j = 0; j < global->ligand->number_of_atoms; j++ ) {
    if(ligand->orbit == ADDED) continue;

    for(i = 0; i < num_targ_atoms; i++){

      sq_dist = global->target_ligand_sq_dists[j*num_targ_atoms + i];
      if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

      if(atoms_overlap2(&target[i], &ligand[j], sq_dist, max_overlap, 
                        &overlap, &distance)){
#ifdef TRACE
        printf("BUMP %2d: %s %s %s (%5.3f) <--> %s %s %s (%5.3f) : %f (%f) "
               "act = %d, %d\n", counter, target->name, target->residue,
               target->residue_num, target->rad, global->ligand_file_name,
               ligand->residue, ligand->residue_num, ligand->rad,
               distance, distance - target->rad - ligand->rad + max_overlap,
               target->act, ligand->act);
#endif
        counter++;
        resolve = FAILURE;
      }
    }
  }

  for ( i = 0; i < global->number_of_waters; i++ ) {
    water = &global->waters[i];
    if ( water->state == CONSERVED ) {
      for ( j = 0; j < global->number_of_target_atoms; j++ ) {
	target = &global->target_atoms[j];
	distance = dist_fun ( water->pos, 
			      target->pos );
	if ( ( target->act == ACCEPTOR
	       || target->act == DONEPTOR
	       || target->act == DONOR )
	     && distance < MIN_HBOND_LENGTH )
	  /* displaced by polar target atom */
	  water->state = POLAR_DISPLACED;
	if ( target->act == NOTHING 
	     && distance < WATER_RAD + target->rad - max_overlap )    
	  /* displaced by non-polar target atom */
	  water->state = DISPLACED;
      }
    }
    if ( water->state == CONSERVED ) {
      for ( j = 0; j < global->ligand->number_of_atoms; j++ ) {
	ligand = &global->ligand->atoms[j];
	if ( ligand->orbit != ADDED ) {
	  distance = dist_fun ( water->pos, ligand->pos );
	  if ( ( ligand->act == ACCEPTOR
		 || ligand->act == DONEPTOR
		 || ligand->act == DONOR )
	       && distance < MIN_HBOND_LENGTH )
	    /* displaced by polar ligand atom */
	    water->state = POLAR_DISPLACED;
	  if ( ligand->act == NOTHING 
	       && distance < WATER_RAD + ligand->rad - max_overlap )  
	    /* displaced by non-polar ligand atom */
	    water->state = DISPLACED;
	}
      }
    }
  }
  return resolve;
}

int  bump_check ( global_data_pt  global )
{
  int i, j, k;
  atom_pt  water;
  bump_t   bump_points[MAX_ANCHOR_BUMPS];
  float distance;
  float sq_dist;
  int number_of_bumps;
  float overlap;
  dist_bin_pt bin = 0;
  molecule_pt ligand = global->ligand;
  atom_pt ligand_atoms = global->ligand->atoms;
  atom_pt ligand_beg = ligand->atoms;
  atom_pt ligand_end = ligand->atoms + ligand->number_of_atoms;
  atom_pt targ_atom, lig_atom;
  atom_pt *targ_atom_p;

  float overlap_tolerance = global->anchor_overlap;
  int *ligand_matches = global->current_orientation.ligand_matches;
  global->number_of_anchor_corrections = 0;

#ifdef TRACE
  printf("\n");
  for(i = 0; i < 3; i++ ) {
    j = ligand_matches[i];
    printf("%3d: %2d %5s %2d %s\n", 
           global->current_orientation.template_matches[i], j,
           ligand_atoms[j].type_str, ligand_atoms[j].fragment,
           global->current_orientation.matched_type[i] == HYDROPHOB ?  "(HPHOB)" : "" );
  }
#endif
  
  for(i = 0; i < 3; i++)
    global->anchor_fragments[i] = ligand_atoms[ligand_matches[i]].fragment;
  find_rigid_bonds(global->ligand, ligand_matches);
  find_anchor_fragments_and_anchor_atoms ( global );

  /* Check for ligand anchor bumps */
  for(number_of_bumps = 1; number_of_bumps > 0; ){
    number_of_bumps = 0;

    for(lig_atom = ligand_beg; lig_atom < ligand_end; ++lig_atom){
      /* toneroma 07MAR07 - do consider only non added H's in bump check */
      /* Makes little sense here since Hs are never ANCHOR atoms */
      if(lig_atom->orbit == ADDED || lig_atom->res != ANCHOR) 
        continue;

      bin = get_bin(&global->target_dists_array, lig_atom->pos);
      if(bin == 0) continue;
      for(targ_atom_p = bin->atoms; targ_atom_p < bin->atoms + bin->num_atoms; 
          ++targ_atom_p){
        targ_atom = *targ_atom_p;
        if(targ_atom->type > 5 && targ_atom->level != HETATM) continue;

        sq_dist = squared_dist(lig_atom->pos, targ_atom->pos);
        if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

        if(atoms_overlap2(targ_atom, lig_atom, sq_dist, overlap_tolerance, 
                          &overlap, &distance)){
          if(number_of_bumps >= MAX_ANCHOR_BUMPS ){
            /* the vector `bump_point' is already full, thus there
               is no hope for this ligand ==> bump-check failed */
            if(BUMP_ANCHOR > global->last ) global->last = BUMP_ANCHOR;
            return FAILURE;
          }

          for(k = 0; k < 3; k++){
            bump_points[number_of_bumps].protein_atom_pos[k] 
              = (double) targ_atom->pos[k];
            bump_points[number_of_bumps].ligand_atom_pos[k] 
              = (double) lig_atom->pos[k];
          }
          bump_points[number_of_bumps].protein_atom_rad =
            (double) targ_atom->rad;
          bump_points[number_of_bumps].ligand_atom_rad =
            (double) lig_atom->rad;
          bump_points[number_of_bumps].overlap = 
            (double) -1.0*(distance - targ_atom->rad - lig_atom->rad);
#ifdef TRACE
          printf("BUMP no. %2d: %3s %s %3s <--> %4s %2d : %5.3f (%5.3f)\n",
                 number_of_bumps+1, targ_atom->name, targ_atom->residue, 
                 targ_atom->residue_num, lig_atom->type_str, 
                 lig_atom->fragment, distance, 
                 targ_atom->rad + lig_atom->rad - overlap_tolerance - 
                 distance);
#endif
          number_of_bumps++;
        }
      }
    }
    if(number_of_bumps > 0){
      if(unbump_anchor(global, bump_points, number_of_bumps) == FAILURE){
#ifdef TRACE
        printf("Unbump anchor failed!\n");
#endif
        if(BUMP_ANCHOR > global->last) global->last = BUMP_ANCHOR;
        return FAILURE;
      }
      global->number_of_anchor_corrections++;
      if(global->number_of_anchor_corrections > MAX_ANCHOR_CORRECTION_STEPS){
#ifdef TRACE
        printf("Tried %d anchor corrections, failed\n", 
               global->number_of_anchor_corrections );
#endif
	if ( BUMP_ANCHOR > global->last ) global->last = BUMP_ANCHOR;
	return FAILURE;
      }
    }
  }

  /* so far we have ignored binding-site waters, now let's see if there is
     any overlap between a water and an atom in the anchor fragment, is so,
     it is displaced */
  for(i = 0; i < global->number_of_waters; i++ ) {
    water = &global->waters[i];
    for(j = 0; j < ligand->number_of_atoms; j++){
      if(ligand->fragment_locations[ligand_atoms[j].fragment] != ANCHOR || 
         ligand_atoms[j].orbit == ADDED) continue;

      distance = dist_fun ( global->waters[i].pos, ligand_atoms[j].pos );
      if(ligand_atoms[j].act == NOTHING && 
         distance < WATER_RAD + ligand_atoms[j].rad - overlap)
        water->state = DISPLACED;
      else if ((ligand_atoms[j].act == ACCEPTOR || 
                ligand_atoms[j].act == DONEPTOR || 
                ligand_atoms[j].act == DONOR)
               && distance < MIN_HBOND_LENGTH) 
        water->state = POLAR_DISPLACED;
    }
  }
  
#ifdef TRACE
  printf("No more anchor bumps after %d correction steps\n", 
	 global->number_of_anchor_corrections );
#endif
#ifdef FILTER_BURIED_CARBONS
  if(check_buried_carbons(global) == FAILURE) return FAILURE;
#endif

  for(i = 0; i < ligand->number_of_flexible_bonds; i++ ) {
    ligand->bond_directions[i] = UNKNOWN;
    if(ligand->bond_types[i] == RIGID){
      j = ligand_atoms[ligand->bonds[ligand->flexible_bonds[i]].atom1].fragment;
      ligand->fragment_locations[j] = ANCHOR;
      j = ligand_atoms[ligand->bonds[ligand->flexible_bonds[i]].atom2].fragment;
      ligand->fragment_locations[j] = ANCHOR;
    }
  }
  
#ifdef TRACE
  for(i = 0; i < ligand->number_of_atoms; i++)
    printf("%2d: %5s %2d frag %2d\n", i + 1, ligand->atoms[i].type_str,
           ligand->atoms[i].number, ligand->atoms[i].fragment );
#endif
  direct_flexible_bonds ( global );
#ifdef TRACE
  ligand_trace(ligand);
#endif

  /* now check, which atoms have to be checked for intramolecular bumps after
     rotation, i.e., identify those atoms with a distance of more than two
     covalent bonds */
  identify_second_grade_neighbors(global);
  if(unbump_side_chains(global) == FAILURE){
    if(global->number_of_mean_field_optimizations > 0 ){
#ifdef TRACE
      printf ( "%d rotations tried in %d iterations, %d bumps left\n",
	       global->number_of_side_chain_rotations,
	       global->number_of_mean_field_optimizations,
	       global->number_of_bumps );
      printf("Total overlap = %5.3f\n", 
             global->current_orientation.total_overlap );
#endif
      if(global->current_orientation.total_overlap < 
         global->finally_tolerated_overlap){
#ifdef TRACE
        print_target_ligand_distances(global->target_atoms, 
                                      global->number_of_target_atoms, 
                                      ligand->atoms, ligand->number_of_atoms);
#endif 
        return SUCCESS;
      }
    }
#ifdef TRACE
    printf("Unresolvable side-chain collision, ligand ruled out\n");
#endif
#ifdef OUTPUT_BUMPS
    find_all_bumps_output ( global );
#endif

    if(BUMP_SIDE_CHAIN > global->last ) global->last = BUMP_SIDE_CHAIN;
    return FAILURE;
  } /*3*/
  
#ifdef TRACE
  printf ( "No more ligand side-chain collisions!\n" );
  printf ( "No target side chain collisions!\n" );
#endif

  /* check if the full bump check is too stringent -- i.e. if something
   * would pass the above method but fail here, we are doing something or
   * someone a disservice */
  if(full_bump_check(global) != SUCCESS ){
#ifdef TRACE
    printf("FINAL BUMP CHECK FAILED for %s_%d\n", global->ligand_file_name,
	   global->binding_modes_counter);
#endif
    return FAILURE;      
  }

  if(global->current_orientation.total_overlap < 
     global->finally_tolerated_overlap){
#ifdef TRACE
    print_target_ligand_distances(global->target_atoms, 
                                  global->number_of_target_atoms, 
                                  ligand->atoms, ligand->number_of_atoms);
#endif 
    return SUCCESS;
  }
  return FAILURE;
}

void ligand_trace(molecule_pt ligand)
{
  int i; 
  atom_pt ligand_atoms = ligand->atoms;

  for ( i = 0; i < ligand->number_of_fragments; i++ )
    if ( ligand->fragment_locations[i] == ANCHOR )
      printf ( "%2d: ANCHOR (dist=%d,way_to_anchor=%d)\n", i,
	       ligand->anchor_dist[i],
	       ligand->way_to_anchor[i] );
    else
      printf("%2d: OFF ANCHOR (dist=%d,way_to_anchor=%d,bond=%d)\n", i,
	     ligand->anchor_dist[i],
	     ligand->way_to_anchor[i],
	     ligand->flexible_bonds[ligand->way_to_anchor[i]]);
  printf ( "flexible bonds: " );
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    printf ( " %d", ligand->flexible_bonds[i] );
  printf ( "\n" );
  printf ( "flexible bonds in anchor: " );
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( ligand->bond_types[i] == RIGID )
      printf ( " %d", ligand->flexible_bonds[i] );
  printf ( "\n" );
  printf ( "Rigid bonds: \n" );
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( ligand->bond_types[i] == RIGID )
      printf("%2d: %2d - %2d\n", ligand->flexible_bonds[i], 
	     ligand->bonds[ligand->flexible_bonds[i]].atom1 + 1,
	     ligand->bonds[ligand->flexible_bonds[i]].atom2 + 1 );
  for ( i = 0; i < ligand->number_of_atoms; i++ )
    printf ( "(%d: %d) ", i + 1, ligand_atoms[i].fragment );
  printf ( "\n" );
  printf ( "Remaining flexible bonds:\n" );
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( ligand->bond_types[i] == FLEX )
      printf("%2d (bond %2d): %2d - %2d %2d %s\n", i,
	     ligand->flexible_bonds[i], 
	     ligand->bonds[ligand->flexible_bonds[i]].atom1 + 1,
	     ligand->bonds[ligand->flexible_bonds[i]].atom2 + 1,
	     ligand->bond_directions[i],
	     ligand->bond_directions[i] == STRAIGHT ?
	     "STRAIGHT" : ligand->bond_directions[i] == REVERSE ?
	     "REVERSE" : ligand->bond_directions[i] == END ?
	     "END" : "UNKNOWN" );
    else
      printf("%2d: %2d - %2d %2d %s (RIGID)\n", 
	     ligand->flexible_bonds[i], 
	     ligand->bonds[ligand->flexible_bonds[i]].atom1 + 1,
	     ligand->bonds[ligand->flexible_bonds[i]].atom2 + 1,
	     ligand->bond_directions[i],
	     ligand->bond_directions[i] == STRAIGHT ?
	     "STRAIGHT" : ligand->bond_directions[i] == REVERSE ?
	     "REVERSE" : ligand->bond_directions[i] == END ?
	     "END" : "UNKNOWN" );
  
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    printf ( "(%d:%d:%d:%d)", i,
	     ligand->flexible_bonds[i],
	     ligand->bond_types[i],
	     ligand->bond_directions[i] );
  printf ( "\n" );
}

void print_target_ligand_distances(atom_pt target, int num_target_atoms,
                                   atom_pt ligand, int num_ligand_atoms)
{
  int i, j;
  float distance;
  for(i = 0; i < num_ligand_atoms; i++ ) {
    for(j = 0; j < num_target_atoms; j++ ){
      distance=dist_fun(ligand[i].pos, target[j].pos);
      printf("L_%d_%s: %9.4f %9.4f %9.4f <-> T_%s_%s: %7.3f %7.3f %7.3f "
             "<=> dist= %f\n", i+1, ligand[i].type_str, ligand[i].pos[0], 
             ligand[i].pos[1], ligand[i].pos[2], target[j].residue, 
             target[j].name, target[j].pos[0], target[j].pos[1], 
             target[j].pos[2], distance);
    }
  }
}
