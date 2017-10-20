#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <types.h>
#include <basics.h>
#include <least_square_fit.h>
#include <dist_fun.h>
#include <err_handle.h>
#include <transform_molecule.h>

#if 0
#endif

int transform_molecule(global_data_pt global, int *template_matches,
                       int *matched_type, int *ligand_matches)
{
  double         rot_matrix[3][3],
    trans_vector[3];
  float err;
  double         new_points[3][3],
    weights[3],
    ligand_points[3][3],
    template_points[3][3];
  interaction_pt template_pt;
  int i,j,k;
  double new_pos[3];

  double rmstol = global->rmstol;
  float *ligand_pos = global->compound_interactions->pos;
  template_pt = global->template_interactions;

  for(i = 0; i < 3; i++ )
    for(j = 0; j < 3; j++ ){
      ligand_points[i][j] = (double) ligand_pos[3*ligand_matches[i]+j];
      template_points[i][j] = (double) template_pt[template_matches[i]].pos[j];
    }

  /* we do a weight superposition, where a hbond point is weight three
     times as important as a hydrophobic interaction center, since
     the matching of the latter needs to be less exact */
  for (i = 0; i < 3; i++)
    if(matched_type[i] == HYDROPHOB ) weights[i] = HPHOB_MATCH_WEIGHT;
    else weights[i] = 1.0;

  if(least_square_fit(ligand_points, template_points, new_points, weights, 
		      rot_matrix, trans_vector ) == FAILURE )
    return NO_MATCH;

#ifdef TRACE
  for(i = 0; i < 3; i++ )
    printf ( "(%5.3lf,%5.3lf,%5.3lf) -> (%5.3lf,%5.3lf,%5.3lf)\n",
	     new_points[i][0],
	     new_points[i][1],
	     new_points[i][2],
	     template_points[i][0],
	     template_points[i][1],
	     template_points[i][2] );
#endif

  err = 0.0;
  /* again, when computing the RMSD for the superposition of the two
     triangles, allow a less exact fit for the hydrophobic points */
  for (i = 0; i < 3; i++)
    err += weights[i]*weights[i] * 
           squared_dist_double(template_points[i], new_points[i]);

  if(err > 3.0 * rmstol*rmstol){
#ifdef TRACE
    printf ( "FAILURE: rmsd = %5.3lf\n", sqrt(err/3.0) );
#endif
    if( RMS_DEVI > global->last ) global->last = RMS_DEVI;
    return NO_MATCH;

  }else{ 
    /* the ligand atoms, the target atoms, and the waters have to be 
       reset to their original states and positions, since they might have
       been rotated or translated for the last mapping */
    memcpy(global->ligand->atom_positions, global->ligand->orig_atom_positions,
           3*global->ligand->number_of_atoms * sizeof(float));
    for(i = 0; i < global->ligand->number_of_atoms; i++ ){
	global->ligand->atoms[i].act = global->ligand->orig_act[i];
	global->ligand->atoms[i].rad = global->ligand->orig_rad[i];
    }
    /* mark all added hydrogens as added, if there was a
       position assigned during scoring, the type_str was
       changed from "H_ADD" to "H" */ 
    for(i = global->ligand->number_of_atoms 
          - global->ligand->number_of_added_hydrogens; 
        i < global->ligand->number_of_atoms; i++ )
      sprintf ( global->ligand->atoms[i].type_str, "H_ADD" );	      

    reset_target_and_waters(global);
  }

  /* if one ligand atom was matched to a template point that has to
     be covalently bonded to the target, its radius has to be reduced
     so that there won't be any trouble during the bump checks */
  for(i = 0; i < 3; i++ )
    if(template_pt[template_matches[i]].key_point == COVALENT )
      global->ligand->atoms[global->compound_interactions->atom_index[ligand_matches[i]]].rad = 0.5; 

#ifdef TRACE
    printf ( "MATCH: rmsd = %5.3lf\n", sqrt(err/3.0) );
#endif
 
    global->rms_err = (float) err;  
    for(i = 0; i < global->ligand->number_of_atoms; i++ ){
      for(j = 0; j < 3; j++ ){
	new_pos[j] = trans_vector[j];
	for(k = 0; k < 3; k++){
	  new_pos[j] += rot_matrix[j][k] * global->ligand->atom_positions[3*i + k];
	}
      }
      for (j = 0; j < 3; j++ )
	global->ligand->atom_positions[3*i + j] = (float) new_pos[j];
    }
  for(i = 0; i < 3; i++ ){
    ligand_matches[i] =
      global->compound_interactions->atom_index[ligand_matches[i]];
      
#ifdef TRACE
    if(matched_type[i] != HYDROPHOB )
      printf("atom %d (%s %s %s %d) matched to template point %d %s\n", 
             ligand_matches[i] + 1,
             global->ligand->atoms[ligand_matches[i]].name,
             global->ligand->atoms[ligand_matches[i]].type_str,
             global->ligand->atoms[ligand_matches[i]].act == DONOR ? "D" : 
             global->ligand->atoms[ligand_matches[i]].act == ACCEPTOR ? "A" : 
             global->ligand->atoms[ligand_matches[i]].act == DONEPTOR ? "AD" : 
             global->ligand->atoms[ligand_matches[i]].act == HYDROPHOB ? "H" :
             "ERROR",
             global->ligand->atoms[ligand_matches[i]].act,
             template_matches[i],
             template_pt[template_matches[i]].act == DONOR ? "D" : 
             template_pt[template_matches[i]].act == ACCEPTOR ? "A" : 
             template_pt[template_matches[i]].act == DONEPTOR ? "AD" : 
             template_pt[template_matches[i]].act == HYDROPHOB ? "H" :
             "ERROR" );
    else
      printf("ring center matched to template point %d %s\n",
             template_matches[i],
             template_pt[template_matches[i]].act == DONOR ? "D" : 
             template_pt[template_matches[i]].act == ACCEPTOR ? "A" : 
             template_pt[template_matches[i]].act == DONEPTOR ? "AD" : 
             template_pt[template_matches[i]].act == HYDROPHOB ? "H" :
             "ERROR" );
#endif

    /* when we matched ligand atoms to template points, the matching was
       based on the number of the atom in the mol2 file, but since there
       might have been some dummy atoms, these numbers have to be updated */
    ligand_matches[i] = global->ligand->atom_index[ligand_matches[i]+1];
  }

  return MATCH;
}


void reset_target_and_waters(global_data_pt global)
{
  int i, j;
  memcpy(global->target_atom_positions, global->orig_target_atom_positions,
         3*global->number_of_target_atoms *
         sizeof(*global->target_atom_positions));
  for ( i = 0; i < global->number_of_target_atoms; i++ )
    global->target_atoms[i].act = global->orig_target_atom_act[i];
  for ( i = 0; i < global->number_of_target_residues; i++ ){
    global->target_rotations[i] = NO;
    global->target_intra_overlap[i] = NO;
  }

  /* reset the waters */
  for ( i = 0; i < global->number_of_waters; i++ ) {
    global->waters[i].state = CONSERVED;
    for ( j = 0; j < 3; j++ )
      global->waters[i].pos[j] = global->orig_water_positions[i][j];
  }
}
