#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <distance_matrices.h>
#include "types.h"
#include "basics.h"
#include "find_all_bumps.h"
#include "compute_all_rotation_angles.h"
#include "compute_unbump_dependencies.h"
#include "initialize_unbump_matrices.h"
#include "rotatable_bonds.h"
#include "trans_rotate.h"
#include "dist_fun.h"
#include "mean_field_minimization.h"
#include "rotate_unbump_bonds.h"
#include "debug_funs.h"

#define FILE_OUTPU

void pre_mean_field_trace(unbump_data_pt **unbump_data, 
                          int number_of_unbump_bonds, int number_of_bumps);

void post_mean_field_trace(unbump_data_pt **unbump_data, int *unbump_bonds,
                           residue_pt target_residues, int **unbump_dependents,
                           int *number_of_unbump_dependents, 
                           int number_of_unbump_bonds, int number_of_bumps);

int  compute_matrix ( atom_pt            side_chain_atoms,
		      atom_pt            ligand_atom,
		      int                res_type,
		      int                bond,
		      transform_matrix_t matrix )
{
  int base; 
  int head;
  
  for(base = 0;
      side_chain_atoms[base].type != rotatable_bonds[res_type][bond-1]; base++);
  for(head = 0; 
      side_chain_atoms[head].type != rotatable_bonds[res_type][bond]; head++);

  /* the rotatable bond is between atoms[base] and atoms[head] */
  compute_transformation_matrix ( ligand_atom->pos,
				  side_chain_atoms[base].pos,
				  side_chain_atoms[head].pos,
				  matrix );
  return side_chain_atoms[head].level; 
}

double compute_unbump_angle(double x2, double d, double y1, double y2,
                            double lxz, double l)
{
  double  angle;
  
  angle = lxz * lxz
    - ( 1.0 / 2.0 / x2
	* ( x2 * x2
	    + y2 * y2
	    - 2 * y1 * y2
	    + l * l
	    - d * d ) )
    * ( 1.0 / 2.0 / x2
	* ( x2 * x2
	    + y2 * y2
	    - 2 * y1 * y2
	    + l * l
	    - d * d ) );
  /* If the angle is negative there is no solution to the equation, 
   * i.e., we cannot compute a rotation angle 
   */
  if(angle < 0.0) return FLT_MAX;
  angle = sqrt ( angle ) / lxz;
  angle = fabs ( asin ( angle ) );
  if ( ( 1.0 / 2.0 / x2
		     * ( x2 * x2
			 + y2 * y2
			 - 2 * y1 * y2
			 + l * l
			 - d * d ) ) < 0 )
    /* new x ordinate of side-chain atom is < 0 */
    angle *= (-1);
  return angle;
}


/*
 *  It might happen that resolving the overlap by rotating around
 *  the current rotation axis is not possible. This usually happens
 *  when the rotated atom is too close to this axis. This routine computes
 *  the maximal obtainable distance when rotating the rotate atom to the
 *  opposite site of the fixed atom.
 */
int  check_resolvability ( float              fixed_atom_pos[3],
			   float              rotated_atom_pos[3],
			   float              final_distance,
			   transform_matrix_t matrix )
{
  float  help_fixed[3],
         help_rotated[3];
  float  angle;
  int    i;

  /* transform the positions of the bumping atoms to the new system
     with target_pos is on the XY-plane */
  for(i = 0; i < 3; i++){
    help_rotated[i] = rotated_atom_pos[i];
    help_fixed[i] = fixed_atom_pos[i];
  }
  transform_point ( help_fixed, matrix );
  transform_point ( help_rotated, matrix );

  /* compute angle between ligand atom position (i.e. the atom which
     will be rotated to resolve the bump) and the X-axis */
  angle =  ( help_rotated[X] * help_rotated[X] ) 
    + ( help_rotated[Z] * help_rotated[Z] );
  /* there has been some trouble with small values in this routine,
     so, whenever it seems like the value is getting too small, just
     assume that there is no hope and leave this routine */
  /* Gotta luv the above comment, would be interesting to know what actually
   * going on rather than just fudge the return if there is "trouble" */
  if ( angle < SMALL_FLOAT ) return BUMP;
  /*  printf("help_rotated %f %f %f\n", help_rotated[X], help_rotated[Y], help_rotated[Z]);*/
  if((help_rotated[X] * help_rotated[X]) < SMALL_FLOAT){
      /*      printf ("small_value\n");*/
      if (help_rotated[Z] < 0 ) angle = -0.5 * M_PI;
      else angle = 0.5 * M_PI;
  }else
    angle = asin(help_rotated[Z] / 
                 sqrt((double) (help_rotated[X] * help_rotated[X] + 
                                help_rotated[Z] * help_rotated[Z])));
  /*  printf("unbump_pre_angle %f\n", angle);*/
  /* compute the angle for getting the maximal distance between both
     points, i.e. rotate 'help' into the XY-plane opposed to the 
     target point */
  if(help_rotated[Z] < 0) angle = (-3) * M_PI - fabs(angle);
  else angle = 3 * M_PI - fabs(angle);
  /*  printf("unbump_angle %f\n", angle);*/
  rotate_single_atom_around_y_axis(help_rotated, angle);

  if(compare_double(squared_dist(help_rotated, help_fixed), 
                    final_distance*final_distance ) == -1 )
    return BUMP;
  return NO_BUMP;
}


int  unbump_side_chains ( global_data_pt  global )
{
  unbump_data_pt **unbump_data;
  int            iteration,
                 result,
                 number_of_actual_rotations;
  int            i;
  char           filename[FILENAME_MAX];

#ifdef FILE_OUTPUT
  char           filename[FILENAME_MAX];
#endif
#ifdef TRACE
  int j;
#endif

/*************************************/

  iteration = 0;
  number_of_actual_rotations = 1;
  result = BUMP;
  /* to make the computation of the intermolecular distances more efficient,
     here we only once compute the distances and store them in a lookup table,
     so that we do not have to call 'dist_fun()' over and over when modeling
     the induced complementarity */
  /* Added the 0.0005 to keep things the same as the previous version */
  initialize_inter_dist_matrix(global->target_atom_positions, 
                               global->number_of_target_atoms,
                               global->ligand->atom_positions,
                               global->ligand->number_of_atoms,
                               global->target_ligand_sq_dists, 
                               DONT_CARE_BUMP_DISTANCE + 0.0005);
  result = find_all_bumps(&global->current_orientation, global);

#ifdef TRACE
  printf("Before unbump: %d bumps (result = %s)\n", 
	 global->number_of_bumps,
	 result == BUMP ? "BUMP" : result == FAILURE ? "FAILURE" :
	 result == NO_BUMP ? "NO_BUMP" : "DON'T KNOW" );
#endif

  global->number_of_mean_field_optimizations = 0;
  global->number_of_side_chain_rotations = 0;
  global->number_of_ligand_side_chain_rotations = 0;
  global->number_of_target_side_chain_rotations = 0;
  global->number_of_water_translations = 0;
  while(iteration < MAX_LIGAND_BUMP_RESOLVE_ITERATIONS && result == BUMP && 
        number_of_actual_rotations > 0){
    global->number_of_unbump_bonds = 0;
    unbump_data = global->unbump_data;

    /* Zero out all the pointers */
    for ( i = 0; i < MAX_SIDE_CHAIN_BUMPS; i++ )
      memset(unbump_data[i], 0, MAX_UNBUMP_BONDS*sizeof(unbump_data[i][0]));

    /* At present compute_all_rotation_angles always returns SUCCESS */
    /* Check for non-resolvable bumps */
    if(compute_all_rotation_angles ( global ) == FAILURE || 
       global->number_of_unbump_bonds == 0 ){
#ifdef TRACE
      printf("number_of_unbump_bonds: %d\n", global->number_of_unbump_bonds );
      printf("leaving unbump_side_chains due to unresolvable bumps\n" );
#endif

#ifdef OUTPUT_BUMPS
      find_all_bumps_output(global);
#endif
/*************************************/
      return FAILURE;
    }

    compute_unbump_dependencies ( global );
    initialize_unbump_matrices ( global );

#ifdef TRACE
    pre_mean_field_trace(global->unbump_data, global->number_of_unbump_bonds,
                         global->number_of_bumps);
#endif

    mean_field_minimization(global);

#ifdef TRACE
    post_mean_field_trace(global->unbump_data, global->unbump_bonds,
                          global->target_residues, global->unbump_dependents,
                          global->number_of_unbump_dependents,
                          global->number_of_unbump_bonds,
                          global->number_of_bumps);
#endif

    /* check if there is at least one bump that cannot be resolved, 
     * if so, leave this routine */
    if(rotate_unbump_bonds(global, &number_of_actual_rotations ) == FAILURE){
#ifdef TRACE 
      printf ( "There is one unresolvable collision\n" );
#endif  
      return FAILURE; 
    }

    global->number_of_bumped_ligands ++ ; 
    global->number_of_side_chain_rotations += number_of_actual_rotations;
    result = find_all_bumps (&global->current_orientation, global );
    iteration++;
    global->number_of_mean_field_optimizations = iteration;

#ifdef TRACE
    printf("after iteration %d: %d bumps (%d rotations)\n", 
           iteration, global->number_of_bumps, number_of_actual_rotations );
#endif

    if((iteration == MAX_LIGAND_BUMP_RESOLVE_ITERATIONS || 
        number_of_actual_rotations == 0 ) && result != NO_BUMP ){
      /* there are still some bumps left and this was the last iteration 
         or there hasn't any rotation occurred */
#ifdef TRACE
      printf ( "max iteration or no rotations in current iteration\n" );
      if(number_of_actual_rotations == 0) printf("No rotations made\n");
      else printf ( "Tried %d iterations\n", iteration );
#endif	  

#ifdef OUTPUT_BUMPS
      find_all_bumps_output ( global );
#endif

#ifdef FILE_OUTPUT
      sprintf ( filename, "%s/%s/%s/%s_ligands/bump_%s_ligand_%04d.mol2",
                global->data_root, global->protein, global->template,
                global->database, global->ligand_file_name, 
		global->number_of_bumped_ligands );
      write_ligand_mol2 ( global, filename );
      sprintf ( filename, "%s/%s/%s/%s_targets/bump_%s_target_%04d.pdb",
                global->data_root, global->protein, global->template,
                global->database, global->ligand_file_name,
		global->number_of_bumped_ligands);
      write_target_pdb ( global, filename );
      global->number_of_bumped_ligands ++;
#endif
/*************************************/
      return FAILURE;
    }
  }

  global->number_of_mean_field_optimizations = iteration;
  return result;
}
 
void
pre_mean_field_trace(unbump_data_pt **unbump_data, int number_of_unbump_bonds,
                     int number_of_bumps)
{
  int i, j;

  printf ( "Unbump angles\n" );
  for(j = 0; j < number_of_unbump_bonds; j++){
    for(i = 0; i < number_of_bumps; i++ )
      if(unbump_data[i][j] != NULL) printf(" %6.3lf", unbump_data[i][j]->angle);
      else printf ( "  *.***" );
      printf ( "\n" );
  }	  
  printf ( "Unbump forces:\n" );
  for(j = 0; j < number_of_unbump_bonds; j++ ){ 
    for(i = 0; i < number_of_bumps; i++ ) 
      if(unbump_data[i][j] != NULL) printf(" %6.3lf", unbump_data[i][j]->force);
      else printf ( "  *.***" ); 
    printf ( "\n" ); 
  }	  
  printf ( "Unbump penalties:\n" ); 
  for ( j = 0; j < number_of_unbump_bonds; j++ ){ 
    for ( i = 0; i < number_of_bumps; i++ ) 
      if(unbump_data[i][j] != NULL) 
        printf(" %6.3lf", unbump_data[i][j]->penalty); 
      else printf ( "  *.***" ); 
    printf ( "\n" ); 
  }	  
  printf ( "Unbump probabilities before mean-field optimization:\n" ); 
  for ( j = 0; j < number_of_unbump_bonds; j++ ){ 
    for ( i = 0; i < number_of_bumps; i++ ) 
      if(unbump_data[i][j] != NULL) 
        printf(" %6.3lf", unbump_data[i][j]->probability); 
      else printf("  *.***" ); 
    printf("\n"); 
  }	
  printf("%d bumps, %d bonds\n", number_of_bumps, number_of_unbump_bonds);
}

void
post_mean_field_trace(unbump_data_pt **unbump_data, int *unbump_bonds,
                      residue_pt target_residues, int **unbump_dependents, 
                      int *number_of_unbump_dependents,
                      int number_of_unbump_bonds, int number_of_bumps)
{
  int i, j;
  printf("Unbump probabilities after mean-field optimization:\n" ); 
  for ( j = 0; j < number_of_unbump_bonds; j++ ){ 
    for(i = 0; i < number_of_bumps; i++ ) 
      if(unbump_data[i][j] != NULL) 
        printf(" %6.3lf", unbump_data[i][j]->probability ); 
      else printf ( "  *.***" ); 
    printf ( "\n" ); 
  }	
  printf ( "Unbump angles\n" ); 
  for ( j = 0; j < number_of_unbump_bonds; j++ ){ 
    for ( i = 0; i < number_of_bumps; i++ ) 
      if ( unbump_data[i][j] != NULL ) 
        printf ( " %6.3lf", unbump_data[i][j]->angle ); 
      else printf ( "  *.***" ); 
    printf ( "\n" ); 
  }	  
  for ( i = 0; i < number_of_unbump_bonds; i++){
    if(unbump_bonds[i] >= 0){      
      printf("bond %2d ", i ); 
      printf("(%s,%d): ", target_residues[unbump_bonds[i] / 10].num, 
             unbump_bonds[i] % 10); 
      for(j = 0; j < number_of_unbump_dependents[i]; j++) 
        printf("%2d ", unbump_dependents[i][j]); 
      printf("\n" ); 
    }else{ 
      printf("bond %2d (%d): ", i, ( (-1) * unbump_bonds[i] ) - 1 ); 
      for(j = 0; j < number_of_unbump_dependents[i]; j++) 
        printf("%2d ", -1*(unbump_bonds[unbump_dependents[i][j]])-1);
      printf ( "\n" ); 
    }
  }
}
