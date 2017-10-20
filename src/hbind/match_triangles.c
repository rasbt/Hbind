#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "mymalloc.h"
#include "match_triangles.h"
#include "basics.h"
#include "hash_class.h"
#include "transform_molecule.h"
#include "bump_check.h"
#include "quicksort.h"
#include "err_handle.h"
#include "dist_fun.h"
#include "read_mol2.h"
#include "sum_charges.h"
#include "analyze_ligand.h"
#include "transform_molecule.h"
#include "docking_features.h"
#include "score_complex.h"
#include "write_ligand_mol2.h"
#include "write_target_pdb.h"
#include "write_waters_pdb.h"

#define AFFINITY_SCORE 62  		/* affinity trained funtion # 66 */
#define ORIENTATION_SCORE 63	/* affinity trained funtion # 66 */


static int permutations[6][3] = { { 0, 1, 2 },
				  { 0, 2, 1 },
				  { 1, 0, 2 },
				  { 1, 2, 0 },
				  { 2, 0, 1 },
				  { 2, 1, 0 } };

static int point_permutations[6][3] = { { 0, 1, 2 },
					{ 1, 0, 2 },
					{ 2, 1, 0 },
					{ 1, 2, 0 },
					{ 2, 0, 1 },
					{ 0, 2, 1 } };

int ligand_is_small(atom_pt atoms, int number_of_atoms);

int load_ligand_atoms(global_data_pt global);

double  dist_fun_array ( float *array,
			 int   index1,
			 int   index2 )
{
  double  dx, dy, dz;

  dx = (double) array[index1+X] - array[index2+X];
  dy = (double) array[index1+Y] - array[index2+Y];
  dz = (double) array[index1+Z] - array[index2+Z];

  /* define each step to be "double" so that we have the similiar
     result and precision on both PC and SGI */

  return (double) sqrt ( (double) dx * dx + (double) dy * dy + (double) dz * dz );
  
}

int  unique ( int  *array, int  size )
{
  int  largest;
  int  i;

  if ( size > 1 ) quicksort_int( array, 0, size - 1 );

  largest = 0;
  for ( i = 1; i < size; i++ )
    if ( array[i] > array[largest] ){
      largest++;
      array[largest] = array[i];
    }
  return largest + 1;
}

int  compute_class ( int  *array, int  *act, int  my_index, int  count )
{
  if ( my_index < 3 ){
    if ( act[my_index] == DONEPTOR){
      act[my_index] = DONOR;
      count = compute_class ( array, act, my_index + 1, count );
      act[my_index] = ACCEPTOR;
      count = compute_class ( array, act, my_index + 1, count );
      act[my_index] = DONEPTOR;
    }else
      count = compute_class ( array, act, my_index + 1, count );  
  }else{
    array[count] = hash_class[act[0]+act[1]+act[2]];
    count++;
  }
  return count;
}


int  check_triangle ( global_data_pt global,
		      triangle_pt    triangle,
		      double         *sides,
		      int            *act )
{
  double  weights[3];
  double  err, val, min_err;
  float   dmetol;
  int     sum, match, min;
  int     i, j;

  dmetol = global->dmetol;
  min = -1;
  min_err = 9999.9;
  /* check all possible matchings */  
  for ( i = 0; i < 6; i++ ){
    match = TRUE;
    for ( j = 0; j < 3 && match == TRUE; j++ ){
      match = FALSE;
      sum = triangle->act[j] + act[point_permutations[i][j]];
      if ( sum == ACCEPTOR_ACCEPTOR 
	   || sum == DONOR_DONOR 
	   || sum == HYDROPHOB_HYDROPHOB 
	   || sum == DONOR_DONEPTOR
	   || sum == ACCEPTOR_DONEPTOR )
	match = TRUE;
    }
    if ( match == TRUE ){
      err = 0.0;
      /* we are doing a weighted least-squares-fit superposition later,
	 so we are already here allowing a looser fit for the DME of
	 the matched triangles, if a ligand center is matched onto
	 a hydrophobic template point */
      weights[0] = weights[1] = weights[2] = 1.0;
      if ( triangle->act[0] == HYDROPHOB ){
	/* triangle dist[0] is the distance between points 0 and 1,
	   dist[1] between 1 and 2, and dist[2] between 2 and 0, so
	   if point 0 is a hydrophobic point, the matching of distances
	   dist[1] and dist[2] has to be less exact, and so on */
	    
	weights[0] = HPHOB_MATCH_WEIGHT;
	weights[2] = HPHOB_MATCH_WEIGHT;
      }
      if ( triangle->act[1] == HYDROPHOB ){
	    
	weights[0] = HPHOB_MATCH_WEIGHT;
	weights[1] = HPHOB_MATCH_WEIGHT;
      }
      if ( triangle->act[2] == HYDROPHOB ){
	weights[1] = HPHOB_MATCH_WEIGHT;
	weights[2] = HPHOB_MATCH_WEIGHT;
      }
	  
      for ( j = 0; j < 3; j++ ){
	val = weights[j] * ( triangle->dist[j] 
			     - sides[permutations[i][j]] );
	err += val * val;
      }
      err = sqrt ( err / 3.0 );
      if ( compare_double( err, min_err ) == -1 ){
	min_err = err;
	min = i;
      }
    }
  }
  if ( min == -1 || compare_double( min_err, dmetol ) ==1 ){
    if ( DME > global->last ) global->last = DME;
    return NO_MATCH;
  }else {
    global->dist_matrix_err = (float) min_err;
    /* return min+1, since NO_MATCH is 0 */
    return min + 1;
  }
}

int compute_triangles(global_data_pt global, 
                      int *array, int my_index,
                      int length, int offset)
{
  double         side[3];
  float          *pos;
  double         perimeter;
  int            actual_index[3];
  int            act[3],
    type_classes[20];
  int            orig, 
    number_of_type_classes,
    max, min, matched_index,
    result;
  int            i, j, c, p, l, s;
  triangle_parameters_t* triangles;
  int template_matches[3];
  int matched_type[MAX_INTER_PAIRS];
  int ligand_matches[3];
  dock_feats_pt features;

  triangles = &(global->triangle_parameters);
  pos = global->compound_interactions->pos;

  /* get the real index in the array containing the interaction points 
     of the molecule */
  for ( i = 0; i < 3; i++ ){
    actual_index[i] = 3 * offset + 3 * array[i];
    act[i] = global->compound_interactions->act[offset+array[i]];
  }
  side[0] = dist_fun_array ( pos, actual_index[0], actual_index[1] );
  side[1] = dist_fun_array ( pos, actual_index[1], actual_index[2] );
  side[2] = dist_fun_array ( pos, actual_index[2], actual_index[0] );

  /* compute perimeter of the corresponding triangle */
  perimeter = side[0] + side[1] + side[2];

#ifdef TRACE
  printf("side = (%.16f, %.16f, %.16f), perimeter = %.16f\n",
         side[0],side[1],side[2],perimeter);
#endif

  /* sort the side lengths of the triangle */
  max = min = 0;
  if ( compare_double( side[1], side[max] ) == 1 ) max = 1;
  else if (compare_double( side[1], side[min] ) == -1 ) min = 1;
  if ( compare_double( side[2], side[max] ) == 1 ) max = 2;
  else if ( compare_double( side[2], side[min] ) == -1 ) min = 2;

  p = UNKNOWN;
  l = UNKNOWN;
  s = UNKNOWN;
  
  /* This got hacked to allow for the "switching" between triangle parameters
   based on the ligand size.  We build the hash table for the template once.
   In addition we could be screening both large and small ligands.  This
   means that the hash table has all the possible triangle sizes and we need
   start from the value triangles->min_long_side_small when looking for a
   bucket.  However, the range is allowed to change based on the size of the
   ligand.*/
  if(perimeter >= triangles->min_perimeter 
     && perimeter < triangles->max_perimeter)
    p = (int)(BUCKETS_PER_A_PERIMETER * (perimeter - triangles->min_perimeter));
  if(side[max] >= triangles->min_long_side
     && side[max] < triangles->max_long_side )
    l = (int)(BUCKETS_PER_A_LENGTH * 
      (side[max] - SMALL_TRIANGLE_MIN_LONGEST_SIDE));
  if(side[min] >= triangles->min_short_side 
     && side[min] < triangles->max_short_side )
    s = (int)(BUCKETS_PER_A_LENGTH * (side[min] - triangles->min_short_side));

#ifdef TRACE
  printf("p = %6d, l = %6d, s = %6d\n", p, l, s);
#endif

  if ( p != UNKNOWN && l != UNKNOWN && s != UNKNOWN ){
    /* match triangle against all corresponding hash-table entries */
    type_classes[0] = hash_class[act[0]+act[1]+act[2]];
    number_of_type_classes = 1;
    if ( type_classes[0] >= 10 ){
      /* there is at least one doneptor point involved */
	
      number_of_type_classes = compute_class ( type_classes, act, 0, 0 );
      number_of_type_classes = unique ( type_classes, 
					number_of_type_classes );
    }
    for ( c = 0; c < number_of_type_classes; c++ )
      for ( i = 0; 
	    i < global->hash_table[type_classes[c]][p][l][s].number_of_entries; 
	    i++ ){

	if(global->hash_table[type_classes[c]][p][l][s].entries[i] == 0){
	  printf("wait for me\n");
	}

	matched_index = check_triangle(global, 
                        global->hash_table[type_classes[c]][p][l][s].entries[i],
                                       side, act );

	if ( matched_index != NO_MATCH ) {
	  matched_index--;
	  for ( j = 0; j < 3; j++ ) {
	    template_matches[j] = 
	      global->hash_table[type_classes[c]][p][l][s].entries[i]->index[j];
	    matched_type[j] = 
	      global->hash_table[type_classes[c]][p][l][s].entries[i]->act[j];
	    ligand_matches[j] = 
	      offset + array[point_permutations[matched_index][j]];
	  }

	  result = transform_molecule(global, template_matches, matched_type,
                                      ligand_matches);
	  if(result == MATCH){
            features = &(global->current_orientation);
            init_features(features);
            memcpy(features->template_matches, template_matches, 
                   3*sizeof(*template_matches));
            memcpy(features->ligand_matches, ligand_matches, 
                   3*sizeof(*ligand_matches));
            memcpy(features->matched_type, matched_type, 
                   MAX_INTER_PAIRS*sizeof(*ligand_matches));
            if(bump_check(global) == SUCCESS) score(features, global, 0);
          }
	  /* there were problems with this ligand file, so skip it */
	  else if ( result == FATAL_FAILURE ) return FATAL_FAILURE;
	}
      }
  }

  if ( my_index > 0 ) {
    orig = array[my_index-1];    /* store original entry at position index-1 */
    while ( array[my_index-1] < array[my_index] - 1 ) {
      /* increase value at position index - 1 as long as it is
	 still smaller than the value at position index */
      array[my_index-1]++;
      /* recursively increase all elements in the array at
	 positions [0..index-2] */
      result = compute_triangles ( global, array, my_index - 1, 
                                   length, offset );
      /* there were problems with this ligand file, so skip it */
      if ( result == FATAL_FAILURE ) return FATAL_FAILURE;
      /***** Added by PCS -- 16-May-01 *****/	  
      else if (global->last < TRIANGLE_MATCH ) global->last = TRIANGLE_MATCH;
      /*************************************/
    }
    /* restore original value at position index-1 */
    array[my_index-1] = orig;
  }
  return SUCCESS;
}

    
void  match_triangles ( global_data_pt  global )
{
  int array[MAX_NUMBER_OF_MOL2_ATOMS];
  int result;
  int i, j;
  char linebuffer[MAX_MOL2_LINELENGTH];
  char *pos;
  int same_conf = 0;
  int rv = 0;
#ifndef OUTPUT_ALL_MATCHES
  char lig_name[FILENAME_MAX + 1];
#endif

  int number_of_compounds = global->compound_interactions->number_of_compounds;
  pts_compound_pt compounds = global->compound_interactions->compounds;
  /* 2006:10:31 */
  /*  printf("\n\t\t *** # of compounds: %d",number_of_compounds);*/

  /* Initialize variables */
  init_features(&global->best_orientation);
  global->last = FATAL_FAILURE;
  global->binding_modes_counter = 0;

  for ( i = 0; i < number_of_compounds; i++ ) {
    strcpy(global->old_ligand_name_noconf, global->ligand->name_noconf);
    global->number_of_screened_compounds++;
    global->compound_index = i;
    if(global->database_num_exist == TRUE){
      /*      fprintf(stderr, "Screened %d molecules, (total screened conformers: "
              "%d => %.3f %%)\r", global->total_num_molecules_noconf, 
              global->number_of_screened_compounds-1, 
              (100*((float)global->number_of_screened_compounds-1) /
	      (float)global->database_num_conformers));*/
      fprintf(stderr, "Screened:%8d(%6.3f%%)||Top:%-7.3f(AVG:%-7.3f,STD:%6.3f) %s\r", global->total_num_molecules_noconf, 
              (100*((float)global->number_of_screened_compounds-1) /
	       (float)global->database_num_conformers),
	      global->best_affiscore_so_far,
	      global->affiscore_mean,
	      global->affiscore_stdd,
	      global->best_affi_name);
    }else{
      /*      fprintf(stderr, "Screened %d molecules, (total screened conformers: %d)"
              "\r", global->total_num_molecules_noconf, 
              global->number_of_screened_compounds-1);*/
      fprintf(stderr, "Screened:%8d||Top:%-7.3f(AVG:%-7.3f,STD:%6.3f) %s\r", 
	      global->total_num_molecules_noconf, 
	      global->best_affiscore_so_far,
	      global->affiscore_mean,
	      global->affiscore_stdd,
	      global->best_affi_name);
    }
    
    strcpy(global->ligand_file_name, compounds[i].name );
    printf("%s\n", global->ligand_file_name );
    fflush(stdout );

    rv = load_ligand_atoms(global);
    if(rv != SUCCESS){
      if(rv != RESTART_SKIP_MOL){
        printf("Reading of %s: failed\n\t-- Skipping ligand."
               "  Check err file.\n", global->ligand_file_name );
        fprintf(stderr, "Reading of %s: failed\n\t-- Skipping ligand."
                "  Check err file.\n", global->ligand_file_name );
      }
      continue;
    }

    /* Can't make a triangle with 2 points */
    if(compounds[i].number_of_points < 3) {
      sprintf(linebuffer, "ligand %s has < 3 interaction points", 
              compounds[i].name );
      err_warning2("match_traingles", linebuffer);
    }	
      
    reset_target_and_waters(global);
    if(ligand_is_small(global->ligand->atoms, global->ligand->number_of_atoms))
      global->triangle_parameters.min_long_side = 
          SMALL_TRIANGLE_MIN_LONGEST_SIDE;
    else 
      global->triangle_parameters.min_long_side = 
          LARGE_TRIANGLE_MIN_LONGEST_SIDE;

    for(j = 0; j < compounds[i].number_of_points; j++ ) array[j] = j;
    for(j = 2, result = SUCCESS; 
        j < compounds[i].number_of_points && result == SUCCESS; j++) {
      array[2] = array[j];
      result = compute_triangles(global, array, 2, 3, compounds[i].first_point);
    }

    /* need to see if next ligand is a conformer of the current ligand */
    same_conf = 0;
    if(global->group_conformers && i < (number_of_compounds - 1) &&
       strncmp(compounds[i+1].name, "singleton", 9)){
      pos = strrchr(compounds[i+1].name, '_');
      if(strncmp(compounds[i].name, compounds[i+1].name, 
                 pos - compounds[i+1].name) == 0)
      same_conf = 1;
    }
      
    if(!same_conf){
      switch(global->last){
      case FATAL_FAILURE:
        printf("%s: failed FATAL FAILURE -- Skipping ligand. "
               "Check err file.\n", global->ligand_file_name );
        break;
      case TRIANGLE_MATCH:
        printf ( "%s: failed triangle match\n", global->ligand_file_name );
        break;
      case DME:
        printf ( "%s: failed DME check\n", global->ligand_file_name );
        break;
      case RMS_DEVI:
        printf ( "%s: failed RMS deviation\n", global->ligand_file_name );
        break;
      case BUMP_ANCHOR:
        printf ( "%s: failed anchor unbumping\n", global->ligand_file_name );
        break;
      case BUMP_SIDE_CHAIN:
        printf("%s: failed side-chain unbumping\n", global->ligand_file_name );
        break;
      case SCORING:
        printf ( "%s: failed scoring\n", global->ligand_file_name );
        break;
      case PASSED:

#ifdef OUTPUT_ALL_MATCHES	
        printf("%s: %d %s\n", global->ligand_file_name,
	       global->binding_modes_counter == 0 ? 
	       1 : global->binding_modes_counter,
	       global->binding_modes_counter == 0 ? 
	       "orientation" : "orientations" );
#else
        if(global->best_orientation.affi_score <= global->score_cutoff){
          if(global->group_conformers)
            sprintf(lig_name, "%s_0000", global->ligand->name_noconf);
          else
            sprintf(lig_name, "%s%s_0000", global->ligand->name_noconf, 
                    global->ligand->conf_number);
          report_docking(&global->best_orientation, 
                         &global->best_orient_positions, lig_name, global);
        }else{
          fprintf(stderr, "Error in logic flow\n");
          /* We shouldn't get here since we should be stuck in the "SCORING"
           * case if we didn't have a good enough docking
           */
        }
#endif
        global->number_of_potential_ligands++;
        break;
      default:
        break;
      }
      /* store the statistics for how far this compound made it */
      global->filter_counter[global->last]++;
      fflush ( stdout );

      /* Reset variables */
      init_features(&global->best_orientation);
      global->last = FATAL_FAILURE;
      global->binding_modes_counter = 0;
    }  
  }
}

int load_ligand_atoms(global_data_pt global)
{
  char filename[FILENAME_MAX];
  char tmpfilename[FILENAME_MAX];
  int len;
  char *pos;
  FILE *MOL2;
  char err_msg[FILENAME_MAX];
  int rv = 0;

  /* Mol2 file */
  if(!strncmp(global->compound_name, "singleton", 9)){
    len = strlen(global->compound_dir) + strlen(global->ligand_file_name) + 7; 
    if(len > FILENAME_MAX){
      err_error2("load_ligand_atoms", "Ligand file name is too long");
      return FATAL_FAILURE;
    }
    sprintf(filename, "%s/%s.mol2", global->compound_dir,
            global->ligand_file_name );
    if((MOL2 = open_mol2(filename)) == NULL) return FATAL_FAILURE; 
    if((rv = read_mol2(MOL2, filename, global, NULL)) != SUCCESS){
      fclose(MOL2);
      return rv; 
    }
    fclose(MOL2);
  
  /* MultiMol2 file */
  }else{
    /* New multi file */
    if(strcmp(global->compound_name, global->old_compound_name)){
      strcpy(global->old_compound_name, global->compound_name);

      /* Get name of multimol2 file -- may this change ? -- seems like this
         part of the naming is an artifact of at most N records per multimol2 */
      strcpy(tmpfilename, global->compound_name);
      pos = strrchr(tmpfilename, '_');
      if(!pos){
        sprintf(err_msg, 
                "Expected to find an '_' in the ligand compound name\n");
        fprintf(stderr, "%s", err_msg);
        err_print(err_msg);
        printf("%s", err_msg);
        return FATAL_FAILURE;
      }
      *pos = 0;
      len = strlen(global->compound_dir) + strlen(tmpfilename) + 7; 
      if(len > FILENAME_MAX){
        err_error2("load_ligand_atoms", "Ligand file name is too long");
        return FATAL_FAILURE;
      }
      sprintf(filename, "%s/%s.mol2", global->compound_dir, tmpfilename);

      if(global->MM2_FILE) fclose(global->MM2_FILE);
      global->MM2_FILE = NULL;
      if((global->MM2_FILE = open_mol2(filename)) == NULL) 
        return FATAL_FAILURE;
    }

    len = strlen(global->compound_dir) + strlen(global->compound_name) + 7; 
    if(len > FILENAME_MAX){
      err_error2("load_ligand_atoms", "Ligand file name is too long");
      return FATAL_FAILURE;
    }
    sprintf(filename, "%s/%s.mol2", global->compound_dir, 
            global->compound_name);
    rv = read_mol2(global->MM2_FILE, filename, global, 
                   global->ligand_file_name);
    if(rv != SUCCESS){
      if(rv == FATAL_FAILURE) fclose(global->MM2_FILE);
      return rv;
    }
  }

  if(analyze_ligand ( global ) == FAILURE){
    sprintf(err_msg, "ERROR: Analysis of ligand structure failed; "
            "skipping ligand %s\n\n", global->ligand_file_name );
    err_print(err_msg);
    fprintf(stderr, "%s", err_msg);
    printf("%s", err_msg);
    return -1;
  }
  sum_charges(global->ligand);

  if(global->ligand_flag != 0) free(global->ligand_flag);
  global->ligand_flag =
    (short *) mymalloc(3 * global->ligand->number_of_atoms * sizeof(short) );

  return SUCCESS;
}

int ligand_is_small(atom_pt atoms, int number_of_atoms)
{
  int i,j;
  for(i = 0; i < number_of_atoms; i++)
    for(j = i+1; j < number_of_atoms; j++)
      if(dist_fun(atoms[i].pos, atoms[j].pos) > SMALL_LIGAND_DIAMETER
         && atoms[i].type != H && atoms[j].type != H  ){ 
        printf("Large ligand\n"); 
        return FALSE;
      }
  printf("Small ligand\n");
  return TRUE;              
}

int score(dock_feats_pt features, global_data_pt global, FILE *good_acts_file)
{
  int i;
  atom_pt atoms; 
  int num_atoms;
  float *pos;
#ifdef OUTPUT_ALL_MATCHES
  char lig_name[FILENAME_MAX + 1];
#else
  int *flag;
#endif     

  /* Score the given complex (orientation */
  score_complex(global, features, good_acts_file);

  /* Keep the dockings with affinity scores meeting the user's threshold.
   * If FILTER_BURIED_CARBONS is defined, keep only those dockings that
   * have at least 50% of the ligand carbon atoms buried in the protein 
   */
  if(features->affi_score > global->score_cutoff 
#ifdef FILTER_BURIED_CARBONS
     || features->buried_carbons < 0.5 
#endif     
#ifndef OUTPUT_ALL_MATCHES
     /* If we are not outputing all matches, keep the best orientation */
     || features->orient_score > global->best_orientation.orient_score
#endif
    ){
    if(SCORING > global->last) global->last = SCORING;
    return SCORING;
  }

  strcpy(features->ligand_name, global->ligand->name);
  strcpy(features->ligand_name_noconf, global->ligand->name_noconf);
  strcpy(features->ligand_conf, global->ligand->conf_number);

  /* copy the bump counts here for now */
  global->last = PASSED;
  features->number_of_bumps = global->number_of_bumps;
  features->num_anchor_corrections =
    global->number_of_anchor_corrections;
  features->num_mean_field_optimizations =
    global->number_of_mean_field_optimizations;
  features->num_side_chain_rotations =
    global->number_of_side_chain_rotations;
  features->num_ligand_side_chain_rotations =
    global->number_of_ligand_side_chain_rotations;
  features->num_target_side_chain_rotations =
    global->number_of_target_side_chain_rotations;

#ifdef OUTPUT_ALL_MATCHES
  sprintf(lig_name, "%s%s_%04d", features->ligand_name_noconf, 
          features->ligand_conf, global->binding_modes_counter);
  report_docking(features, NULL, lig_name, global);
  global->binding_modes_counter++;
#else
  global->binding_modes_counter++;

  
  /* NOTE: if the dock_feats_t structure is changed to contain pointers,
   * a copy routine will need to be implemented so that a deep copy is
   * performed */
  global->best_orientation = *features;

  /* Copy ligand atom positions */
  atoms = global->ligand->atoms;
  num_atoms = global->ligand->number_of_atoms - 
    global->ligand->number_of_added_hydrogens;
  pos = global->best_orient_positions.ligand_positions;
  for(i = 0; i < num_atoms; i++, pos += 3)
    memcpy(pos, atoms[i].pos, 3*sizeof(*pos)); 

  /* Copy target atom positions */
  atoms = global->target_atoms;
  num_atoms = global->number_of_target_atoms;
  flag = global->best_orient_positions.target_rotations;
  memcpy(flag, global->target_rotations, 
         global->number_of_target_residues*sizeof(*flag)); 
  pos = global->best_orient_positions.target_positions;
  for(i = 0; i < num_atoms; i++, pos += 3)
    memcpy(pos, atoms[i].pos, 3*sizeof(*pos));

  /* Copy water positions + status */
  atoms = global->waters;
  num_atoms = global->number_of_waters;
  flag = global->best_orient_positions.water_states;
  pos = global->best_orient_positions.water_positions;
  for(i = 0; i < num_atoms; ++i, pos += 3){
    memcpy(pos, atoms[i].pos, 3*sizeof(*pos));
    flag[i] = atoms[i].state;
  }

#endif
  return PASSED;
}

void report_docking(dock_feats_pt features, moved_positions_pt moved_positions,
                    char* lig_name, global_data_pt global)
{
  FILE *fp;
  char path_prefix[FILENAME_MAX + 1];
  char mol_fname[FILENAME_MAX + 1];
  char target_fname[FILENAME_MAX + 1];
  char waters_fname[FILENAME_MAX + 1]; 
  char filename[FILENAME_MAX + 1];

  float temp;
  int   change_marker;
  int   errsv;
  char  err_msg[2*FILENAME_MAX];
  
  change_marker = 0;

  write_features_line(features, lig_name, stdout, DOCK_AND_SCORE);
  
  global->total_num_output_molecules++;

  if (features->affi_score < global->best_affiscore_so_far)
    {
      global->best_affiscore_so_far = features->affi_score;
      strcpy(global->best_affi_name, global->ligand->name_noconf);
      change_marker = 1;
    }


  temp = global->affiscore_mean;
  global->affiscore_mean += (features->affi_score - global->affiscore_mean) / (global->total_num_output_molecules);
  global->affiscore_svar += (features->affi_score - temp) * (features->affi_score - global->affiscore_mean);
  global->affiscore_stdd = sqrt(global->affiscore_svar/(global->total_num_output_molecules - 1));

  global->affiscore_significance = (global->best_affiscore_so_far - global->affiscore_mean) / global->affiscore_stdd;
  
  sprintf(filename, "%s/%s/%s/log/%s.status", global->data_root, global->protein, global->template, global->database);
  
  
  fp = fopen (filename, "w");
  
  if (fp == NULL) {
    errsv = errno;
    sprintf(err_msg, "Could not open the file: %s\n\t%s", filename,
	    strerror(errsv));
    fprintf(stderr, err_msg);
    err_print(err_msg);
    err_panic2 ( "write_status_file", "file open failed");
  }
  fprintf(fp, "Number of output molecules: %d, best affiscore: %0.3f :: sig: %0.3f (AVG: %0.3f, STD: %0.3f)\n", global->total_num_output_molecules, global->best_affiscore_so_far, global->affiscore_significance, global->affiscore_mean, global->affiscore_stdd);
  
  fclose (fp);
  
  sprintf(path_prefix, "%s/%s/%s", global->data_root, global->protein,
          global->template);
  sprintf(mol_fname, "%s/%s_ligands/%s.mol2", path_prefix, global->database,
          lig_name);
  sprintf(target_fname, "%s/%s_targets/%s.pdb", path_prefix, global->database,
          lig_name);
  sprintf(waters_fname, "%s/%s_waters/%s.pdb", path_prefix, global->database,
          lig_name);

  if(moved_positions){
    write_ligand_mol2(mol_fname, moved_positions->ligand_positions, features, 
                      global);
    write_target_pdb(global->target_residues, global->number_of_target_residues,
                     global->target_atoms, global->number_of_target_atoms,
                     target_fname, moved_positions->target_rotations,
                     moved_positions->target_positions);
    write_waters_pdb(global->waters, global->number_of_waters, waters_fname, 
                     moved_positions->water_states,  
                     moved_positions->water_positions);
  }else{
    write_ligand_mol2(mol_fname, NULL, features, global);
    write_target_pdb(global->target_residues, global->number_of_target_residues,
                     global->target_atoms, global->number_of_target_atoms,
                     target_fname, global->target_rotations, NULL);
    write_waters_pdb(global->waters, global->number_of_waters, waters_fname, 
                     NULL, NULL);
  }
}
